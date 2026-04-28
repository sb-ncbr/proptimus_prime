# Builder image
FROM python:3.12-slim AS builder

WORKDIR /opt/proptimus

# Setup OS
RUN python3 -m venv /opt/venv

# activate virtual environment
ENV PATH="/opt/venv/bin:$PATH"

# copy files with python requirements and nodejs packages
COPY requirements.txt .

# install python libs
RUN pip install -r requirements.txt

# Runtime image
FROM python:3.12-slim AS base

WORKDIR /opt/proptimus

# Setup OS
RUN apt-get update \ 
    && apt-get install -y --no-install-recommends procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --shell /bin/bash user \
    && chown -R user:user /opt/proptimus

# Switch to non-root user
USER user

# Get prepared python environment
COPY --from=builder /opt/venv /opt/venv

ENV PATH="/opt/venv/bin:$PATH"

# Copy application
COPY prime.py executor_prime.py data_prime.py .

CMD ["python3", "prime.py"]
