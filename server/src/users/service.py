from typing import Union
from uuid import UUID

from sqlalchemy import and_, select, update
from sqlalchemy.ext.asyncio import AsyncSession
from src.users.models import User


class UserService:
    def __init__(self, db_session: AsyncSession):
        self.db_session = db_session

    async def create_user(
        self, username: str, email: str, hashed_password: str, is_admin: bool = False
    ) -> User:
        new_user = User(
            username=username,
            email=email,
            hashed_password=hashed_password,
            is_admin=is_admin,
        )
        self.db_session.add(new_user)
        await self.db_session.flush()
        return new_user

    async def delete_user(self, user_id: UUID) -> Union[User, None]:
        query = (
            update(User)
            .where(and_(User.user_id == user_id, User.is_active is True))
            .values(is_ative=True)
            .returning(User.user_id)
        )
        result = await self.db_session.execute(query)
        deleted_user = result.fetchone()
        return deleted_user[0] if deleted_user else None

    async def get_user(self, user_id: UUID) -> Union[User, None]:
        query = select(User).where(User.user_id == user_id)
        result = await self.db_session.execute(query)
        user = result.fetchone()
        return user[0] if user else None
