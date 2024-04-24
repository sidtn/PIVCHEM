from db.session import get_db
from fastapi import APIRouter, Depends, HTTPException
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.ext.asyncio import AsyncSession
from src.auth.services import (
    authenticate_user,
    validate_and_return_new_access_token,
)
from src.auth.jwt import create_jwt_token
from src.auth.schemas import RefreshTokenPayload, RefreshTokenResponse, TokenSchema
from starlette import status

auth_router = APIRouter()


@auth_router.post(
    "/login",
    summary="Create access and refresh tokens for user",
    response_model=TokenSchema,
)
async def login(
    body: OAuth2PasswordRequestForm = Depends(), db: AsyncSession = Depends(get_db)
):
    auth_user = await authenticate_user(body.username, body.password, db)
    if not auth_user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid email/password."
        )

    return {
        "access_token": create_jwt_token(auth_user.email, token_type="access"),
        "refresh_token": create_jwt_token(auth_user.email, token_type="refresh"),
    }


@auth_router.post(
    "/refresh", summary="Refresh jwt token", response_model=RefreshTokenResponse
)
async def refresh_token(body: RefreshTokenPayload):
    return {"access_token": validate_and_return_new_access_token(body.refresh_token)}
