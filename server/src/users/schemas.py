from pydantic import BaseModel, EmailStr


class ShowUser(BaseModel):
    user_id: int
    email: EmailStr
    is_active: bool


class UserCreate(BaseModel):
    email: EmailStr
    password: str


class DeleteUserResponse(BaseModel):
    deleted_user_id: int


class UpdatedUserResponse(BaseModel):
    updated_user_id: int
