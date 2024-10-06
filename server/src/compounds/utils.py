from io import BytesIO
from tempfile import NamedTemporaryFile
from typing import Union

from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.Draw import MolToImage
from sqlalchemy import select, update
from sqlalchemy.ext.asyncio import AsyncSession
from src.compounds.models import Compound
from src.core.minio import MinioClientFactory

minio_client = MinioClientFactory().get_client()


async def get_next_reg_number(db_session: AsyncSession) -> str:
    query = select(Compound.reg_number).order_by(Compound.reg_number.desc())
    result = await db_session.execute(query)
    last_reg_number = result.first()
    if not last_reg_number:
        return "iboch0000001"
    next_num = int(last_reg_number[0].lstrip("iboch")) + 1
    return "iboch" + str(next_num)


def get_mol_from_smile(smile: str) -> Mol:
    return Chem.MolFromSmiles(smile)


async def generate_img_url_from_mol(mol: Mol, db: AsyncSession) -> None:
    img = BytesIO()
    MolToImage(mol).save(img, format="PNG")
    img_name = f"{hash(Chem.MolToSmiles(mol))}.png"
    res = minio_client.put_image_to_bucket(img_name, BytesIO(img.getvalue()))
    await db.execute(
        update(Compound)
            .where(Compound.structure == mol)
            .values(image_url=res)
        )
