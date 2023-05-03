import uuid

from db.base_class import Base
from sqlalchemy import UUID, Boolean, Column, String


class User(Base):
    user_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    email = Column(String, nullable=False, unique=True)
    is_active = Column(Boolean, default=True)
    is_admin = Column(Boolean, default=False)
    hashed_password = Column(String, nullable=False)

    def __repr__(self):
        return self.email
