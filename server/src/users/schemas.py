import uuid

from pydantic import BaseModel, EmailStr


class ShowUser(BaseModel):
    user_id: uuid.UUID
    email: EmailStr
    is_active: bool


class UserCreate(BaseModel):
    email: EmailStr
    password: str


class DeleteUserResponse(BaseModel):
    deleted_user_id: uuid.UUID


class UpdatedUserResponse(BaseModel):
    updated_user_id: uuid.UUID
