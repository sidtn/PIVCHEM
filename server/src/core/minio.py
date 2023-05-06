from io import BytesIO

from minio import Minio
from settings import MINIO_ACCESS_KEY, MINIO_ENDPOINT, MINIO_SECRET_KEY


class MinioClient:
    def __init__(
        self, endpoint: str, access_key: str, secret_key: str, bucket_name: str
    ):
        self.client = Minio(
            endpoint=endpoint,
            access_key=access_key,
            secret_key=secret_key,
            secure=False,
        )
        self.bucket_name = bucket_name
        if not self.client.bucket_exists(bucket_name):
            self.client.make_bucket(bucket_name)

    def put_image_to_bucket(self, img_name: str, img: BytesIO) -> str:
        self.client.put_object(
            self.bucket_name,
            img_name,
            img,
            length=-1,
            content_type="image/png",
            part_size=5 * 1024 * 1024,
        )
        return self.client.get_presigned_url("GET", self.bucket_name, img_name)


class MinioClientFactory:
    @classmethod
    def get_client(cls):
        return MinioClient(
            endpoint=MINIO_ENDPOINT,
            access_key=MINIO_ACCESS_KEY,
            secret_key=MINIO_SECRET_KEY,
            bucket_name="img-bucket",
        )
