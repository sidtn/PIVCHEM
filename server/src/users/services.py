from typing import Union

from fastapi import HTTPException
from sqlalchemy import and_, select, update
from sqlalchemy.ext.asyncio import AsyncSession
from starlette import status

from src.auth.hasher import Hasher
from src.users.models import User
from src.users.schemas import UserCreate, ShowUser


class UserService:
    def __init__(self, db_session: AsyncSession):
        self.db_session = db_session

    async def get_user_by_email(self, email: str) -> Union[User, None]:
        query = select(User).where(User.email == email)
        result = await self.db_session.execute(query)
        user = result.fetchone()
        return user[0] if user else None

    async def create_user(self, body: UserCreate) -> ShowUser:
        new_user = User(
            email=body.email,
            hashed_password=Hasher.get_password_hash(body.password),
            is_admin=body.is_admin,
        )

        self.db_session.add(new_user)
        await self.db_session.flush()

        return ShowUser(
            user_id=new_user.user_id,
            email=new_user.email,
            is_admin=new_user.is_admin,
            is_active=new_user.is_active
        )

    async def delete_user(self, user_id: int) -> Union[User, None]:
        print(user_id)
        query = (
            update(User)
            .where(and_(User.is_active.is_(True), User.user_id == user_id))
            .values(is_active=False)
            .returning(User.user_id)
        )
        result = await self.db_session.execute(query)
        deleted_user = result.fetchone()
        return deleted_user[0] if deleted_user else None
