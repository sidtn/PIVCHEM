from envparse import Env

env = Env()


DATABASE_URL: str = env.str(
    "REAL_DATABASE_URL",
    default="postgresql+asyncpg://postgres:postgres@0.0.0.0:5432/postgres",
)

SECRET_KEY: str = env.str("SECRET_KEY", default="secret_key")
ACCESS_TOKEN_EXPIRE_MINUTES: int = env.int("ACCESS_TOKEN_EXPIRE_MINUTES", default=30)
