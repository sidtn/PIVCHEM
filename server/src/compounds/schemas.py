from datetime import datetime
from uuid import UUID

from pydantic import BaseModel


class UploadSDFResponse(BaseModel):
    id: int
    smile: str
    reg_number: str
    img_url: str


class CompoundIn(BaseModel):
    reg_number: str
    name: str
    structure: str
    formula: str
    mol_weight: float
    tags: str
    mnr: str
    ms: str
    hplc: str
    mp: str
    doi: str
    cas: str
    storage: str
    img_url: str
    in_stock: bool


class CompoundOut(CompoundIn):
    id: UUID
    created_at: datetime
