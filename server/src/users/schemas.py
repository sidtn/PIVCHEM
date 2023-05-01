import uuid

from pydantic import BaseModel, EmailStr


class ShowUser(BaseModel):
    user_id: uuid.UUID
    username: str
    email: EmailStr
    is_active: bool


class UserCreate(BaseModel):
    username: str
    email: EmailStr
    password: str


class DeleteUserResponse(BaseModel):
    deleted_user_id: uuid.UUID


class UpdatedUserResponse(BaseModel):
    updated_user_id: uuid.UUID
