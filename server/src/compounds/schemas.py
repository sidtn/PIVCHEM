from uuid import UUID

from pydantic import BaseModel


class UploadSDFResponse(BaseModel):
    id: int
    smile: str
    img_url: str


class Compound(BaseModel):
    id: UUID
    reg_number: str
    name: str
    structure: str
    storage: str
