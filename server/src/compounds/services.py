from io import BytesIO
from tempfile import NamedTemporaryFile
from typing import Union

from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.Draw import MolToImage
from src.core.minio import MinioClientFactory

minio_client = MinioClientFactory().get_client()


def generate_img_url_from_mol(mol: Mol) -> Union[str, None]:
    img = BytesIO()
    MolToImage(mol).save(img, format="PNG")
    img_name = f"{hash(Chem.MolToSmiles(mol))}.png"
    res = minio_client.put_image_to_bucket(img_name, BytesIO(img.getvalue()))
    return res


def get_data_from_sdf(sdf_data: NamedTemporaryFile) -> list:
    with Chem.ForwardSDMolSupplier(sdf_data) as sdf:
        mols = [mol for mol in sdf if mol]

    return [
        {
            "id": num,
            "smile": Chem.MolToSmiles(mol),
            "img_url": generate_img_url_from_mol(mol),
        }
        for num, mol in enumerate(mols, start=1)
    ]
