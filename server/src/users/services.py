from datetime import datetime, timedelta
from typing import Union, Sequence

from fastapi import HTTPException
from jose import jwt, ExpiredSignatureError, JWTError
from sqlalchemy import select, update
from sqlalchemy.ext.asyncio import AsyncSession

from settings import SECRET_KEY, ALGORITHM, INVITE_CODE_EXPIRE_DAYS
from src.auth.exceptions import token_expired_exception, credentials_exception
from src.auth.hasher import Hasher
from src.users.models import User
from src.users.schemas import CreateUser, UpdateUser, ActivateUser


class UserService:
    def __init__(self, db_session: AsyncSession):
        self.db_session = db_session

    async def get_user_by_email(self, email: str) -> Union[User, None]:
        query = select(User).where(User.email == email)
        result = await self.db_session.execute(query)
        user = result.fetchone()
        return user[0] if user else None

    @staticmethod
    def _generate_invite_code(user_email: str) -> str:
        expires_at = datetime.now() + timedelta(days=INVITE_CODE_EXPIRE_DAYS)
        to_encode = {"exp": expires_at, "user_email": user_email}
        encoded_jwt = jwt.encode(to_encode, SECRET_KEY, ALGORITHM)
        return encoded_jwt

    async def create_user(self, body: CreateUser) -> User:
        invite_code = self._generate_invite_code(body.email)

        new_user = User(
            email=body.email,
            is_admin=body.is_admin,
            invite_code=invite_code
        )

        self.db_session.add(new_user)
        await self.db_session.flush()

        return new_user

    async def update_user(self, user_id, body: UpdateUser) -> Union[User, None]:
        query = (
            update(User)
            .where(User.user_id == user_id)
            .values(**body.dict())
            .returning(User)
        )
        result = await self.db_session.execute(query)
        user = result.fetchone()
        return user[0] if user else None

    async def activate_user(self, body: ActivateUser) -> Union[User, None]:
        try:
            key = jwt.decode(body.invite_key, SECRET_KEY, algorithms=[ALGORITHM])
            user = await self.get_user_by_email(key.get("user_email"))
            if user.is_active:
                raise HTTPException(status_code=400, detail="Already activated")
        except ExpiredSignatureError:
            raise token_expired_exception
        except JWTError:
            raise credentials_exception

        hashed_password = Hasher.get_password_hash(body.password1)

        query = (
            update(User)
            .where(User.email == key["user_email"])
            .values(hashed_password=hashed_password, is_active=True)
            .returning(User)
        )
        result = await self.db_session.execute(query)
        user = result.fetchone()
        return user[0]

    async def get_users(self) -> Sequence[User]:
        result = await self.db_session.execute(select(User))
        users = result.scalars().all()
        return users
