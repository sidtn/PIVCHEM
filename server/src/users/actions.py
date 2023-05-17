from typing import Union

from src.auth.hasher import Hasher
from src.users.schemas import ShowUser, UserCreate
from src.users.service import UserService


async def _create_new_user(body: UserCreate, session) -> ShowUser:
    async with session.begin():
        user_service = UserService(session)
        user = await user_service.create_user(
            email=body.email,
            hashed_password=Hasher.get_password_hash(body.password),
        )
        return ShowUser(
            user_id=user.user_id,
            email=user.email,
            is_active=user.is_active,
        )


async def _delete_user(user_id, session) -> Union[int, None]:
    async with session.begin():
        user_service = UserService(session)
        deleted_user_id = await user_service.delete_user(
            user_id=user_id,
        )
        return deleted_user_id
