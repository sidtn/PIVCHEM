import uvicorn
from fastapi import APIRouter, FastAPI
from src.auth.handlers import auth_router
from src.compounds.handlers import compounds_router

app = FastAPI(title="Pivchem")

main_api_router = APIRouter()
main_api_router.include_router(auth_router, prefix="/auth", tags=["auth"])
main_api_router.include_router(
    compounds_router, prefix="/compounds", tags=["compounds"]
)
app.include_router(main_api_router)


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8080)
