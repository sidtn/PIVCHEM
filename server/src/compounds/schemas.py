from datetime import datetime

from pydantic import BaseModel, field_serializer
from rdkit import Chem


class UploadSDFResponse(BaseModel):
    smile: str
    props: dict | None


class CompoundIn(BaseModel):
    name: str | None
    structure: str
    storage: str | None = None
    in_stock: bool = False

    @classmethod
    def from_rdkit(cls, compound):
        return cls(
            name=compound.name,
            structure=Chem.MolToSmiles(compound.structure),
            props=compound.props,
            storage=compound.storage,
            in_stock=compound.in_stock,
        )

    def to_rdkit(self):
        return Chem.MolFromSmiles(self.structure)


class CompoundOut(CompoundIn):
    compound_id: int
    reg_number: str | None = None
    created_at: datetime
    image_url: str | None
    props: dict | None = None
    added_by: int | None = None
