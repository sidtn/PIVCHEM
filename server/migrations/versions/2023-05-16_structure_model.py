"""structure model

Revision ID: a9eacc7fc131
Revises: 5882f773b170
Create Date: 2023-05-16 23:04:29.285716

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'a9eacc7fc131'
down_revision = '5882f773b170'
branch_labels = None
depends_on = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.add_column('compound', sa.Column('formula', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('mol_weight', sa.Float(), nullable=True))
    op.add_column('compound', sa.Column('tags', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('mnr', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('ms', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('hplc', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('mp', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('doi', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('cas', sa.String(), nullable=True))
    op.add_column('compound', sa.Column('in_stock', sa.Boolean(), nullable=True))
    # ### end Alembic commands ###


def downgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_column('compound', 'in_stock')
    op.drop_column('compound', 'cas')
    op.drop_column('compound', 'doi')
    op.drop_column('compound', 'mp')
    op.drop_column('compound', 'hplc')
    op.drop_column('compound', 'ms')
    op.drop_column('compound', 'mnr')
    op.drop_column('compound', 'tags')
    op.drop_column('compound', 'mol_weight')
    op.drop_column('compound', 'formula')
    # ### end Alembic commands ###