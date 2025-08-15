import adss
import getpass

from PIL import Image
from astropy.io import fits
import io

class SplusdataError(Exception):
    """Custom exception for S-PLUS data errors."""
    pass


def open_image(image_bytes):
    """Open image from bytes and return PIL Image."""
    from PIL import Image
    im = Image.open(io.BytesIO(image_bytes))
    return im

def save_image(image_bytes, filename):
    """Save image bytes to file."""
    im = open_image(image_bytes)
    im.save(filename)
    

# field frame
class Core:
    """Core interface for interacting with S-PLUS data via ADSSClient."""
    
    def __init__(self, username=None, password=None, SERVER_IP=f"https://splus.cloud", auto_renew=False):
        """Initialize Core with optional credentials and server IP."""
        if not username:
            username = input("splus.cloud username: ")
        if not password:    
            password = getpass("splus.cloud password: ")
            if self.auto_renew:
                self.password = password
        
        self.client = adss.ADSSClient(
            SERVER_IP,
            username=username,
            password=password
        )
        self.collections = []
        
    def _load_collections(self):
        """Load image collections from server."""
        collections = self.client.get_image_collections()
        self.collections = collections

    def get_collection_id_by_pattern(self, pattern):
        """Return first collection whose name contains given pattern."""
        self._load_collections()
        for col in self.collections:
            if pattern in col['name']:
                return col
        raise SplusdataError("Collection not found")
    
    def field_frame(self, field, band, weight=False, outfile=None, _data_release="dr4"):
        """Retrieve a full field frame FITS image for a given field and band."""
        collection = self.get_collection_id_by_pattern(_data_release)
        collection_id = collection['id']
        
        candidates = self.client.list_image_files(collection_id, object_name=field, filter_name=band)

        if len(candidates) == 0 and ("-" in field or "_" in field):
            field = field.replace("-", "_") if "-" in field else field.replace("_", "-")
            candidates = self.client.list_image_files(collection_id, object_name=field)
        if len(candidates) == 0:
            raise SplusdataError(f"Field {field} not found in band {band}")
        
        patterns = collection['patterns']

        final_candidate = None
        f_candidates = []
        
        if weight:
            pattern = "weight"
        for c in candidates:
            if pattern in c['filename']:
                f_candidates.append(c)

        if len(f_candidates) == 0:
            final_candidate = candidates[0]
        elif len(f_candidates) == 1:
            final_candidate = f_candidates[0]
        else:
            fz_candidates = [c for c in f_candidates if c['file_type'] == "fz"]
            final_candidate = fz_candidates[0] if fz_candidates else f_candidates

        image_bytes = self.client.download_image(
            final_candidate['id'],
            output_path=outfile
        )
        
        return fits.open(io.BytesIO(image_bytes))
                                   
    def stamp(self, ra, dec, size, band, weight=False, field_name=None, size_unit="pixels", outfile=None, _data_release="dr4"):
        """Retrieve a FITS cutout (stamp) at given coordinates or field name."""
        collection = self.get_collection_id_by_pattern(_data_release)
        collection_id = collection['id']
        
        if weight:
            weight = "weight"
        if not field_name:
            stamp_bytes = self.client.create_stamp_by_coordinates(
                collection_id=collection_id,
                filter=band,
                ra=ra,
                dec=dec,
                size=size,
                size_unit=size_unit,
                pattern=weight if weight else "",
                output_path=outfile
            )
        else:
            stamp_bytes = self.client.stamp_images.create_stamp_by_object(
                collection_id=collection_id,
                object_name=field_name,
                filter_name=band,
                ra=ra,
                dec=dec,
                size=size,
                size_unit=size_unit,
                pattern=weight if weight else "",
                output_path=outfile
            )
        
        return fits.open(io.BytesIO(stamp_bytes))

    def lupton_rgb(self, ra, dec, size, R="I", G="R", B="G", Q=8, stretch=3, field_name=None, size_unit="pixels", outfile=None, _data_release="dr4"):
        """Retrieve a Lupton RGB composite image."""
        collection = self.get_collection_id_by_pattern(_data_release)
        collection_id = collection['id']

        if not field_name:
            stamp_bytes = self.client.create_rgb_image_by_coordinates(
                collection_id=collection_id,
                ra=ra,
                dec=dec,
                size=size,
                r_filter=R,
                g_filter=G,
                b_filter=B,
                Q=Q,
                size_unit=size_unit,
                stretch=stretch,
                output_path=outfile
            )
        else:
            stamp_bytes = self.client.lupton_images.create_rgb_by_object(
                collection_id=collection_id,
                object_name=field_name,
                ra=ra,
                dec=dec,
                size=size,
                r_filter=R,
                g_filter=G,
                b_filter=B,
                Q=Q,
                size_unit=size_unit,
                stretch=stretch,
                output_path=outfile
            )

        return Image.open(io.BytesIO(stamp_bytes))

    def trilogy_image(self, ra, dec, size, R=["R", "I", "F861", "Z"], G=["G", "F515", "F660"], B=["U", "F378", "F395", "F410", "F430"], noiselum=0.15, satpercent=0.15, colorsatfac=2, size_unit="pixels", field_name=None, outfile=None, _data_release="dr4"):
        """Retrieve a Trilogy RGB composite image."""
        collection = self.get_collection_id_by_pattern(_data_release)
        collection_id = collection['id']

        if not field_name:
            stamp_bytes = self.client.trilogy_images.create_trilogy_rgb_by_coordinates(
                collection_id=collection_id,
                ra=ra,
                dec=dec,
                size=size,
                r_filters=R,
                g_filters=G,
                b_filters=B,
                size_unit=size_unit,
                noiselum=noiselum,
                satpercent=satpercent,
                colorsatfac=colorsatfac,
                output_path=outfile
            )
        else:
            stamp_bytes = self.client.trilogy_images.create_trilogy_rgb_by_object(
                collection_id=collection_id,
                object_name=field_name,
                ra=ra,
                dec=dec,
                size=size,
                r_filters=R,
                g_filters=G,
                b_filters=B,
                noiselum=noiselum,
                size_unit=size_unit,
                satpercent=satpercent,
                colorsatfac=colorsatfac,
                output_path=outfile
            )

        return Image.open(io.BytesIO(stamp_bytes))
    
    def query(self, query, table_upload=None, table_name=None):
        """Run a query on the server, optionally uploading a table."""
        if table_upload is None and table_name is None:
            import pandas as pd
            from astropy.table import Table
            
            table_upload_bytes = None
            if isinstance(table_upload, pd.DataFrame):
                table_upload_bytes = table_upload.to_csv(index=False).encode()
            elif isinstance(table_upload, Table):
                table_upload_bytes = table_upload.to_pandas().to_csv(index=False).encode()
            else:
                raise ValueError("table_upload must be a pandas DataFrame or an astropy Table")

        response = self.client.query_and_wait(
            query_text=query,
            table_name=table_name,
            file=table_upload_bytes
        )
        return response.data