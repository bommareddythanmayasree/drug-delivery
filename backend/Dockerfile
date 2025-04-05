# Use miniconda base image to simplify RDKit install
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy your app files in# Use official Python base with conda (required for RDKit)
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment.yml first for caching
COPY environment.yml .

# Create conda environment with RDKit and dependencies
RUN conda env create -f environment.yml

# Make the environment activate by default
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Copy all files (after env is built)
COPY . .

# Expose port Render expects
EXPOSE 8080

# Run unzip and start the Flask app using gunicorn
CMD ["conda", "run", "--no-capture-output", "-n", "myenv", "bash", "-c", "unzip -o model.zip && gunicorn app:app --bind 0.0.0.0:8080"]
to the container
COPY . /app

# Install RDKit using conda
RUN conda install -y -c rdkit rdkit \
    && pip install --upgrade pip \
    && pip install -r requirements.txt

# Expose port Render expects
EXPOSE 8080

# Start your app using gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "app:app"]
