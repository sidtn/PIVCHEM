import json
from io import BytesIO

from minio import Minio
from settings import MINIO_ACCESS_KEY, MINIO_ENDPOINT, MINIO_SECRET_KEY


class MinioClient:
    def __init__(
        self, endpoint: str, access_key: str, secret_key: str, bucket_name: str
    ):
        self.endpoint = endpoint
        self.bucket_name = bucket_name
        self.client = Minio(
            endpoint=self.endpoint,
            access_key=access_key,
            secret_key=secret_key,
            secure=False,
        )
        self.bucket_name = bucket_name
        if not self.client.bucket_exists(bucket_name):
            self.client.make_bucket(bucket_name)
            policy = {
                "Version": "2012-10-17",
                "Statement": [
                    {
                        "Effect": "Allow",
                        "Principal": {"AWS": "*"},
                        "Action": ["s3:GetBucketLocation", "s3:ListBucket"],
                        "Resource": f"arn:aws:s3:::{bucket_name}",
                    },
                    {
                        "Effect": "Allow",
                        "Principal": {"AWS": "*"},
                        "Action": "s3:GetObject",
                        "Resource": f"arn:aws:s3:::{bucket_name}*",
                    },
                ],
            }
            self.client.set_bucket_policy(bucket_name, json.dumps(policy))

    def _get_img_url(self, object_name):
        url = f"http://{self.endpoint}/{self.bucket_name}/{object_name}"
        return url

    def put_image_to_bucket(self, img_name: str, img: BytesIO) -> str:
        self.client.put_object(
            self.bucket_name,
            img_name,
            img,
            length=-1,
            content_type="image/png",
            part_size=5 * 1024 * 1024,
        )
        return self._get_img_url(img_name)


class MinioClientFactory:
    @classmethod
    def get_client(cls):
        return MinioClient(
            endpoint=MINIO_ENDPOINT,
            access_key=MINIO_ACCESS_KEY,
            secret_key=MINIO_SECRET_KEY,
            bucket_name="img-bucket",
        )
