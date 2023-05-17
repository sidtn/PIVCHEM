from db.base_class import Base
from sqlalchemy import Boolean, Column, Integer, String


class User(Base):
    user_id = Column(Integer, primary_key=True)
    email = Column(String, nullable=False, unique=True)
    is_active = Column(Boolean, default=True)
    is_admin = Column(Boolean, default=False)
    hashed_password = Column(String, nullable=False)

    def __repr__(self):
        return self.email
