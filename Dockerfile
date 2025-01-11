FROM bielcardona/freecad:1.0.0

RUN pip install poetry==1.8.4

WORKDIR /otsun

COPY pyproject.toml poetry.lock README.md ./
COPY src ./src


RUN poetry config virtualenvs.create false
RUN poetry install

