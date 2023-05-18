from typing import Annotated, List

from db.session import get_db
from fastapi import APIRouter, Depends, HTTPException, UploadFile
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from src.auth.actions import user_is_admin
from src.compounds.models import Compound
from src.compounds.schemas import CompoundIn, CompoundOut, UploadSDFResponse
from src.compounds.service import CompoundService
from src.compounds.utils import get_data_from_sdf
from src.users.models import User
from starlette import status


compounds_router = APIRouter(dependencies=[Depends(user_is_admin)])


@compounds_router.post(
    "/upload_sdf",
    summary="Upload sdf file",
    response_model=List[UploadSDFResponse],
)
async def upload_sdf(file: UploadFile):
    if not file.filename.endswith(".sdf"):
        raise HTTPException(
            status_code=status.HTTP_406_NOT_ACCEPTABLE,
            detail="File extension must be sdf",
        )
    try:
        result = get_data_from_sdf(file.file)
        return result
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


# @compounds_router.get("/")
# async def get_compounds(db: AsyncSession = Depends(get_db)):
#     query = select(Compound)
#     result = await db.execute(query)
#     structure = result.all()
#     return {"structure": structure}
