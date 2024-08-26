from datetime import datetime, timedelta
from typing import Any, Union

from jose import jwt
from settings import (
    ACCESS_TOKEN_EXPIRE_MINUTES,
    ALGORITHM,
    REFRESH_TOKEN_EXPIRE_MINUTES,
    SECRET_KEY,
)


def create_jwt_token(user_id: int, subject: str,  token_type: str = "access") -> str:
    delta = (
        ACCESS_TOKEN_EXPIRE_MINUTES
        if token_type == "access"
        else REFRESH_TOKEN_EXPIRE_MINUTES
    )
    expires_at = datetime.now() + timedelta(minutes=delta)

    to_encode = {"exp": expires_at, "id": user_id, "sub": subject, "token_type": token_type}
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, ALGORITHM)
    return encoded_jwt
