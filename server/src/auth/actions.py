from sqlalchemy.ext.asyncio import AsyncSession
from src.users.service import UserService


async def _get_user_by_email(email: str, session: AsyncSession):
    async with session.begin():
        user_service = UserService(session)
        return await user_service.get_user_by_email(
            email=email,
        )
