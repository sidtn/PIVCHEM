from typing import Union

from db.session import get_db
from fastapi import Depends, HTTPException
from fastapi.security import OAuth2PasswordBearer
from jose import ExpiredSignatureError, JWTError, jwt
from settings import ALGORITHM, SECRET_KEY
from sqlalchemy.ext.asyncio import AsyncSession
from src.auth.exceptions import (
    credentials_exception,
    invalid_refresh_token,
    token_expired_exception,
)
from src.auth.hasher import Hasher
from src.auth.jwt import create_jwt_token
from src.users.models import User
from src.users.services import UserService
from starlette import status

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/auth/signin")


async def get_user_by_email(email: str, session: AsyncSession) -> Union[User, None]:
    async with session.begin():
        user_service = UserService(session)
        return await user_service.get_user_by_email(
            email=email,
        )


async def authenticate_user(
    email: str, password: str, db: AsyncSession
) -> Union[User, None]:
    user = await get_user_by_email(email=email, session=db)
    if not user:
        return
    if not Hasher.verify_password(password, user.hashed_password):
        return
    return user


async def get_user_by_jwt_token(
    token: str = Depends(oauth2_scheme),
    db: AsyncSession = Depends(get_db),
) -> Union[User, None]:
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        email = payload.get("sub")
        token_type = payload.get("token_type")
        if not email or token_type != "access":
            raise credentials_exception
    except ExpiredSignatureError:
        raise token_expired_exception
    except JWTError:
        raise credentials_exception

    user = await get_user_by_email(email, db)
    if not user or not user.is_active:
        raise credentials_exception
    return user


def user_is_admin(user: User = Depends(get_user_by_jwt_token)) -> Union[User, None]:
    if not user.is_admin:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Forbidden for non admin users",
        )
    return user


def validate_and_return_new_access_token(token: str) -> Union[str, None]:
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        email = payload.get("sub")
        token_type = payload.get("token_type")
        if not email or token_type != "refresh":
            raise invalid_refresh_token
    except ExpiredSignatureError:
        raise token_expired_exception
    except JWTError:
        raise invalid_refresh_token

    return create_jwt_token(email, token_type="access")
