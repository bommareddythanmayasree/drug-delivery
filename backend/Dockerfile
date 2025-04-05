# Use Miniconda base image for easier RDKit installation
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment file first (for caching)
COPY environment.yml .

# Create conda environment with RDKit and other dependencies
RUN conda env create -f environment.yml

# Activate conda environment in subsequent RUN/CMD steps
SHELL ["conda", "run", "--no-capture-output", "-n", "myenv", "/bin/bash", "-c"]

# Copy all files after environment is built
COPY . .

# Install extra Python packages from requirements.txt (if needed)
RUN pip install --upgrade pip && pip install -r requirements.txt

# Expose port used by Render
EXPOSE 8080

# Unzip the model and run the app with Gunicorn
CMD ["conda", "run", "--no-capture-output", "-n", "myenv", "bash", "-c", "unzip -o model.zip && gunicorn app:app --bind 0.0.0.0:8080"]
