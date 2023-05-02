from fastapi import HTTPException
from starlette import status

credentials_exception = HTTPException(
    status_code=status.HTTP_401_UNAUTHORIZED,
    detail="Incorrect authentication credentials.",
)
token_expired_exception = HTTPException(
    status_code=status.HTTP_401_UNAUTHORIZED, detail="Token is expired"
)
invalid_refresh_token = HTTPException(
    status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid refresh token"
)
