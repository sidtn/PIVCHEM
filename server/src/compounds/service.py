from fastapi.encoders import jsonable_encoder
from sqlalchemy.ext.asyncio import AsyncSession
from src.compounds.models import Compound
from src.compounds.schemas import CompoundIn, CompoundOut
from src.compounds.utils import (
    generate_img_url_from_mol,
    get_mol_from_smile,
    get_next_reg_number,
)
from src.users.models import User


class CompoundService:
    def __init__(self, db_session: AsyncSession):
        self.db_session = db_session

    async def create_compound(self, obj_in: CompoundIn, added_by: User) -> CompoundOut:
        obj_in_data = jsonable_encoder(obj_in)
        smile = obj_in_data.get("structure")
        mol = get_mol_from_smile(smile)
        next_reg_number = await get_next_reg_number(self.db_session)
        obj_in_data["reg_number"] = next_reg_number
        obj_in_data["image_url"] = generate_img_url_from_mol(mol)

        new_compound = Compound(**obj_in_data, added_by=added_by)
        self.db_session.add(new_compound)
        await self.db_session.flush()

        return CompoundOut(
            **obj_in_data,
            id=new_compound.compound_id,
            added_by=added_by.user_id,
            created_at=new_compound.created_at
        )
