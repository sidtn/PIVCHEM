from fastapi import Depends, APIRouter, HTTPException
from sqlalchemy.ext.asyncio import AsyncSession
from starlette import status

from db.session import get_db
from src.auth.services import user_is_admin
from src.users.schemas import ShowUser, UserCreate
from src.users.services import UserService

user_router = APIRouter(dependencies=[Depends(user_is_admin)])


@user_router.post("/", summary="User creation", response_model=ShowUser)
async def create_user(body: UserCreate, db: AsyncSession = Depends(get_db)) -> ShowUser:
    user_service = UserService(db)
    if await user_service.get_user_by_email(body.email):
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=f"user {body.email} already exists",
        )
    return await user_service.create_user(body)


@user_router.delete("/{user_id}", summary="User removing")
async def delete_user(user_id: int, db: AsyncSession = Depends(get_db)):
    user_service = UserService(db)
    deleted_user = await user_service.delete_user(user_id)
    if not deleted_user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="not found",
        )
    return {"user_id": user_id}
