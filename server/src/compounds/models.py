from datetime import datetime

from db.base_class import Base

# from razi.rdkit_postgresql.types import Mol
from sqlalchemy import Boolean, Column, DateTime, Float, ForeignKey, Integer, String
from sqlalchemy.orm import relationship


class Compound(Base):
    compound_id = Column(Integer, primary_key=True)
    reg_number = Column(String, nullable=False, unique=True)
    name = Column(String)
    structure = Column(String)
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

    # __table_args__ = (
    #     Index("compounds_structure", "structure", postgresql_using="gist"),
    # )

    def __repr__(self):
        return f"{self.email} - {self.structure}"
