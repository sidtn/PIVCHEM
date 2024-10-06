from datetime import datetime

from sqlalchemy.dialects.postgresql import JSONB

from db.base_class import Base
from rdkit import Chem
from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    func,
    TypeDecorator,
)
from sqlalchemy.orm import relationship


class Mol(TypeDecorator):
    impl = String

    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        if isinstance(value, Chem.Mol):
            return Chem.MolToSmiles(value)
        raise ValueError(f"Unexpected value: {value}")

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        return Chem.MolFromSmiles(value)


class Compound(Base):
    compound_id = Column(Integer, primary_key=True)
    reg_number = Column(String, nullable=True, unique=True)
    name = Column(String, nullable=True)
    structure = Column(Mol, nullable=False, unique=True)
    mol_weight = Column(Float, nullable=True)
    mnr = Column(String, nullable=True)
    ms = Column(String, nullable=True)
    hplc = Column(String, nullable=True)
    mp = Column(String, nullable=True)
    doi = Column(String, nullable=True)
    cas = Column(String, nullable=True)
    storage = Column(String, nullable=True)
    added_by_id = Column(Integer, ForeignKey("users.user_id"))
    added_by = relationship("User", backref="compounds")
    created_at = Column(DateTime, default=datetime.now)
    image_url = Column(String, nullable=True)
    props = Column(JSONB, nullable=True)
    in_stock = Column(Boolean, default=False)

    __table_args__ = (
        Index("compounds_structure", "structure", postgresql_using="gist"),
    )

    def __repr__(self):
        return f"{self.structure}"
