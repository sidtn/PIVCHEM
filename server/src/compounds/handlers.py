import urllib
from typing import Annotated

from rdkit import Chem
from sqlalchemy import select

from db.session import get_db
from fastapi import APIRouter, Depends, HTTPException, UploadFile
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.ext.asyncio import AsyncSession
from src.auth.services import user_is_admin
from src.compounds.models import Compound
from src.compounds.schemas import CompoundIn, CompoundOut
from src.compounds.services import CompoundService
from src.users.models import User
from starlette import status


compounds_router = APIRouter(dependencies=[Depends(user_is_admin)])


@compounds_router.post(
    "/upload_sdf",
    summary="Upload sdf file",
    response_model=None,
)
async def upload_sdf(
        file: UploadFile,
        user: Annotated[User, Depends(user_is_admin)],
        db: AsyncSession = Depends(get_db),
):
    if not file.filename.endswith(".sdf"):
        raise HTTPException(
            status_code=status.HTTP_406_NOT_ACCEPTABLE,
            detail="File extension must be sdf",
        )
    try:
        inserted_count = 0
        chunk_size = 5000
        mols = []
        with Chem.ForwardSDMolSupplier(file.file) as sdf:
            for mol in sdf:
                if len(mols) < chunk_size:
                    mols.append(mol)
                else:
                    async with db.begin():
                        data = [
                            {
                                "structure": mol,
                                "props": mol.GetPropsAsDict(),
                                "added_by_id": user.user_id,
                            } for mol in mols if mol is not None
                        ]
                        stmt = insert(Compound).values(data).on_conflict_do_nothing()
                        await db.execute(stmt)
                        mols = []
                        inserted_count += len(data)
        return {'detail': {'compounds_added': inserted_count}}
    except Exception as err:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"File decoding error {err}"
        )


@compounds_router.post("/", response_model=CompoundOut)
async def add_compound(
    compound: CompoundIn,
    user: Annotated[User, Depends(user_is_admin)],
    db: AsyncSession = Depends(get_db),
) -> CompoundOut:
    new_compound = await CompoundService(db).create_compound(compound, user)
    return new_compound


@compounds_router.get("/")
async def get_compound(contains: str, db: AsyncSession = Depends(get_db)) -> CompoundOut:
    encoded_smiles = urllib.parse.quote(contains)
    print(encoded_smiles)
    mol = Chem.MolFromSmiles(encoded_smiles)
    print(mol)
    query = select(Compound).filter(Compound.structure == mol)
    result = await db.execute(query)
    compound = result.scalars().all()
    if not compound:
        raise HTTPException(status_code=404, detail="Item not found")
    print(compound[0].__dict__)
    return {}



# @compounds_router.get("/")
# async def get_compounds(contains: str, limit: int = 10, db: AsyncSession = Depends(get_db)) -> dict[
#     str, list[CompoundOut]]:
#     query = select(Compound).filter(Compound.structure == contains)
#     result = await db.execute(query)
#     return {"results": [CompoundOut(**res.__dict__) for res in result.scalars().all()]}

