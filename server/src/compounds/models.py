import uuid
from datetime import datetime

from db.base_class import Base
from razi.rdkit_postgresql.types import Mol
from sqlalchemy import UUID, Column, DateTime, ForeignKey, Index, String
from sqlalchemy.orm import relationship


class Compound(Base):
    compound_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    reg_number = Column(String, nullable=False, unique=True)
    name = Column(String)
    structure = Column(Mol, nullable=False)
    storage = Column(String, nullable=False)
    added_by_uuid = Column(UUID(as_uuid=True), ForeignKey("user.user_id"))
    added_by = relationship("User", back_populates="compounds")
    created_at = Column(DateTime, default=datetime.utcnow)
    image_url = Column(String)

    __table_args__ = (
        Index("compounds_structure", "structure", postgresql_using="gist"),
    )

    def __repr__(self):
        return f"{self.email} - {self.structure}"
