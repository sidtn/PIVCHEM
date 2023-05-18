from datetime import datetime

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
)
from sqlalchemy.orm import relationship
from sqlalchemy.sql.type_api import UserDefinedType


class Mol(UserDefinedType):
    cache_ok = True

    def get_col_spec(self, **kw):
        return "mol"

    def bind_processor(self, dialect):
        def process(value):
            if isinstance(value, Chem.Mol):
                value = Chem.MolToSmiles(value)
            return value
        return process

    def column_expression(self, col):
        return func.mol_to_pkl(col, type_=self)

    def result_processor(self, dialect, coltype):
        def process(value):
            return Chem.MolToSmiles(Chem.Mol(value)) if value else None
        return process


class Compound(Base):
    compound_id = Column(Integer, primary_key=True)
    reg_number = Column(String, nullable=False, unique=True)
    name = Column(String)
    structure = Column(Mol, nullable=False, unique=True)
    formula = Column(String)
    mol_weight = Column(Float)
    tags = Column(String)
    mnr = Column(String)
    ms = Column(String)
    hplc = Column(String)
    mp = Column(String)
    doi = Column(String)
    cas = Column(String)
    storage = Column(String, nullable=False)
    added_by_id = Column(Integer, ForeignKey("user.user_id"))
    added_by = relationship("User", backref="compounds")
    created_at = Column(DateTime, default=datetime.utcnow)
    image_url = Column(String)
    in_stock = Column(Boolean, default=True)

    __table_args__ = (
        Index("compounds_structure", "structure", postgresql_using="gist"),
    )

    def __repr__(self):
        return f"{self.reg_number}"
