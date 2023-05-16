from typing import List

from fastapi import APIRouter, Depends, HTTPException, UploadFile
from src.auth.actions import user_is_admin
from src.compounds.schemas import UploadSDFResponse, CompoundIn, CompoundOut
from src.compounds.services import get_data_from_sdf
from starlette import status

compounds_router = APIRouter(dependencies=[Depends(user_is_admin)])


@compounds_router.post(
    "/upload_sdf", summary="Upload sdf file",
    response_model=List[UploadSDFResponse]
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


@compounds_router.post("/",  response_model=CompoundOut)
async def add_compound(body: CompoundIn):
    return body
