from pydantic import BaseModel


class TokenSchema(BaseModel):
    access_token: str
    refresh_token: str


class RefreshTokenPayload(BaseModel):
    refresh_token: str


class RefreshTokenResponse(BaseModel):
    access_token: str
