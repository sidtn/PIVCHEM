from envparse import Env

env = Env()


DATABASE_URL: str = env.str(
    "REAL_DATABASE_URL",
    default="postgresql+asyncpg://postgres:postgres@0.0.0.0:5432/postgres",
)

SECRET_KEY: str = env.str("SECRET_KEY", default="secret_key")
ACCESS_TOKEN_EXPIRE_MINUTES: int = env.int("ACCESS_TOKEN_EXPIRE_MINUTES", default=1000)
REFRESH_TOKEN_EXPIRE_MINUTES: int = env.int(
    "REFRESH_TOKEN_EXPIRE_MINUTES", default=60 * 24 * 7
)
ALGORITHM: int = env.str("ALGORITHM", default="HS256")

MINIO_ENDPOINT: str = env.str("MINIO_ENDPOINT", default="localhost:9000")
MINIO_ACCESS_KEY: str = env.str("MINIO_ACCESS_KEY", default="admin")
MINIO_SECRET_KEY: str = env.str("MINIO_SECRET_KEY", default="minio_secret_key")
