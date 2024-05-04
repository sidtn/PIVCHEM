from pydantic import BaseModel, EmailStr, model_validator, field_validator
from typing_extensions import Self


class ShowUser(BaseModel):
    user_id: int
    email: EmailStr
    is_admin: bool
    is_active: bool


class CreateUser(BaseModel):
    email: EmailStr
    is_admin: bool = False


class UpdateUser(BaseModel):
    is_admin: bool
    is_active: bool


class ActivateUser(BaseModel):
    invite_key: str
    password1: str
    password2: str

    @field_validator("password1")
    @classmethod
    def check_password_complexity(cls, password: str) -> str:
        errors = ""
        if len(password) < 8:
            errors += "Password must be at least 8 characters long. "
        if not any(character.islower() for character in password):
            errors += "Password should contain at least one lowercase character. "
        if not any(character.isupper() for character in password):
            errors += "Password should contain at least one uppercase character."
        if not any(character.isdigit() for character in password):
            errors += "Password should contain at least one digit."
        if errors:
            raise ValueError(errors)
        return password

    @model_validator(mode="after")
    def check_passwords_match(self) -> Self:
        pw1 = self.password1
        pw2 = self.password2
        if pw1 != pw2:
            raise ValueError("Passwords do not match")
        return self
