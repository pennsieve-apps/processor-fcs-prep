FROM --platform=linux/amd64 python:3.11.16-slim-bookworm

WORKDIR /processor

# Pinned Python dependencies (binary wheels; FlowIO is pure Python from its sdist; no apt packages).
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir --prefer-binary -r /processor/requirements.txt \
    && find /usr/local/lib/python3.11 -depth -type d -name __pycache__ -exec rm -rf {} + \
    && rm -rf /usr/local/lib/python3.11/site-packages/pyarrow/include /usr/local/lib/python3.11/site-packages/pyarrow/tests

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    HOME=/tmp

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.main"]
