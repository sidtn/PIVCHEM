import uvicorn
from fastapi import APIRouter, FastAPI
from src.auth.handlers import auth_router

app = FastAPI(title="Pivchem")

main_api_router = APIRouter()
main_api_router.include_router(auth_router, prefix="/auth", tags=["auth"])
app.include_router(main_api_router)


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8080)
