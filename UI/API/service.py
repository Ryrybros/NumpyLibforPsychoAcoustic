from fastapi import FastAPI, Body
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

import json
from Models.model import model


app = FastAPI()
message = "Hello World"
import os
print(f"currently working at : {os.getcwd()}")

origins = [
    "http://localhost",
    "http://localhost:8080",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/home")
async def index():
    return {"message": message}


class RequestBody(BaseModel):
    title: str
    signal : list[float]
    fs : int
    model : str
    free: bool



@app.post("/home/sendAudio")
async def create_item(value: RequestBody = Body(...)):
    global message
    message = value.title
    mod = model(value.free)
    if(value.model == "moore1997"):
        return dict(mod.moore1997(earSig=value.signal, fs = value.fs))
    else:
        res = mod.glasberg2002(inSig=value.signal, fs = value.fs)
        LTL = res.LTL
        STL = res.STL
        return {
            "LTL": LTL,
            "STL" : STL
        }

