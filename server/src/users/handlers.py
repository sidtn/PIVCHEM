from fastapi import Depends, APIRouter, HTTPException, BackgroundTasks
from sqlalchemy.ext.asyncio import AsyncSession
from starlette import status

from db.session import get_db
from src.auth.jwt import create_jwt_token
from src.auth.services import user_is_admin
from src.users.schemas import ShowUser, CreateUser, UpdateUser, ActivateUser
from src.users.services import UserService

user_router = APIRouter()


@user_router.post(
    "/", summary="User creation", dependencies=[Depends(user_is_admin)], response_model=ShowUser, status_code=201)
async def create_user(
        body: CreateUser, background_tasks: BackgroundTasks, db: AsyncSession = Depends(get_db)) -> ShowUser:
    user_service = UserService(db)
    if await user_service.get_user_by_email(body.email):
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"user {body.email} already exists",
        )
    user = await user_service.create_user(body)

    background_tasks.add_task(print, user.invite_code)

    return ShowUser(
        user_id=user.user_id,
        email=user.email,
        is_admin=user.is_admin,
        is_active=user.is_active
    )


@user_router.patch(
    "/{user_id}",
    summary="User removing",
    dependencies=[Depends(user_is_admin)],
    response_model=UpdateUser,
    status_code=202
)
async def update_user(user_id: int, body: UpdateUser, db: AsyncSession = Depends(get_db)):
    user_service = UserService(db)
    user = await user_service.update_user(user_id, body)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="not found",
        )
    return UpdateUser(
        is_admin=user.is_admin,
        is_active=user.is_active
    )


@user_router.put("/activate", summary="Activate user", status_code=202)
async def activate_user(body: ActivateUser, db: AsyncSession = Depends(get_db)):
    user_service = UserService(db)
    user = await user_service.activate_user(body)
    return {
        "access_token": create_jwt_token(user.email, token_type="access"),
        "refresh_token": create_jwt_token(user.email, token_type="refresh"),
    }
