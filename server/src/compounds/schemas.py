from datetime import datetime

from pydantic import BaseModel


class UploadSDFResponse(BaseModel):
    id: int
    smile: str
    reg_number: str
    image_url: str


class CompoundIn(BaseModel):
    name: str
    structure: str
    formula: str
    tags: str
    storage: str
    in_stock: bool


class CompoundOut(CompoundIn):
    id: int
    reg_number: str
    created_at: datetime
    image_url: str
    added_by: int
