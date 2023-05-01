from database import get_db
from fastapi import APIRouter, Depends, HTTPException
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.exc import IntegrityError
from sqlalchemy.ext.asyncio import AsyncSession
from src.auth.actions import _get_user_by_email
from src.auth.hasher import Hasher
from src.auth.jwt import create_jwt_token
from src.auth.schemas import TokenSchema
from src.users.actions import _create_new_user
from src.users.schemas import ShowUser, UserCreate
from starlette import status

auth_router = APIRouter()


@auth_router.post("/signup", summary="User registration", response_model=ShowUser)
async def create_user(body: UserCreate, db: AsyncSession = Depends(get_db)) -> ShowUser:
    try:
        if await _get_user_by_email(body.email, db):
            raise HTTPException(
                status_code=status.HTTP_409_CONFLICT,
                detail={"detail": f"user {body.email} already exists"},
            )
        return await _create_new_user(body, db)
    except IntegrityError as err:
        raise HTTPException(status_code=503, detail=f"Server error: {err}")


@auth_router.post(
    "/signin",
    summary="Create access and refresh tokens for user",
    response_model=TokenSchema,
)
async def login(
    body: OAuth2PasswordRequestForm = Depends(), db: AsyncSession = Depends(get_db)
):
    user = await _get_user_by_email(body.username, db)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail={"detail": "Incorrect email or password"},
        )
    verified = Hasher.verify_password(body.password, user.hashed_password)

    if not verified:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail={"detail": "Incorrect email or password"},
        )

    return {
        "access_token": create_jwt_token(user.email, token_type="access"),
        "refresh_token": create_jwt_token(user.email, token_type="refresh"),
    }
