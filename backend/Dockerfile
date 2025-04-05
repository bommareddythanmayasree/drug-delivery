# Use miniconda base image to simplify RDKit install
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy your app files into the container
COPY . /app

# Install RDKit using conda
RUN conda install -y -c rdkit rdkit \
    && pip install --upgrade pip \
    && pip install -r requirements.txt

# Expose port Render expects
EXPOSE 8080

# Start your app using gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "app:app"]
