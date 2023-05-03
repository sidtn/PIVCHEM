"""initial

Revision ID: 5882f773b170
Revises: 
Create Date: 2023-05-03 23:52:07.493072

"""
import razi
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '5882f773b170'
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('user',
    sa.Column('user_id', sa.UUID(), nullable=False),
    sa.Column('email', sa.String(), nullable=False),
    sa.Column('is_active', sa.Boolean(), nullable=True),
    sa.Column('is_admin', sa.Boolean(), nullable=True),
    sa.Column('hashed_password', sa.String(), nullable=False),
    sa.PrimaryKeyConstraint('user_id'),
    sa.UniqueConstraint('email')
    )

    op.execute("CREATE EXTENSION rdkit")

    op.create_table('compound',
    sa.Column('compound_id', sa.UUID(), nullable=False),
    sa.Column('reg_number', sa.String(), nullable=False),
    sa.Column('name', sa.String(), nullable=True),
    sa.Column('structure', razi.rdkit_postgresql.types.Mol(), nullable=False),
    sa.Column('storage', sa.String(), nullable=False),
    sa.Column('added_by_uuid', sa.UUID(), nullable=True),
    sa.Column('created_at', sa.DateTime(), nullable=True),
    sa.Column('image_url', sa.String(), nullable=True),
    sa.ForeignKeyConstraint(['added_by_uuid'], ['user.user_id'], ),
    sa.PrimaryKeyConstraint('compound_id'),
    sa.UniqueConstraint('reg_number')
    )
    op.create_index('compounds_structure', 'compound', ['structure'], unique=False, postgresql_using='gist')
    # ### end Alembic commands ###


def downgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_index('compounds_structure', table_name='compound', postgresql_using='gist')
    op.drop_table('compound')
    op.drop_table('user')
    # ### end Alembic commands ###