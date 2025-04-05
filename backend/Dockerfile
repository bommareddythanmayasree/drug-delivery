v# Use Miniconda base image for easier RDKit installation
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Install unzip for extracting model.zip
RUN apt-get update && apt-get install -y unzip

# Copy environment file first (for caching)
COPY environment.yml .

# Create conda environment with RDKit and dependencies
RUN conda env create -f environment.yml

# Activate conda env for following RUN/CMD
SHELL ["conda", "run", "--no-capture-output", "-n", "myenv", "/bin/bash", "-c"]

# Copy all source code after env is built
COPY . .

# Install Python packages
RUN pip install --upgrade pip && pip install -r requirements.txt

# Expose port used by Render
EXPOSE 8080

# Unzip model and start Flask app with Gunicorn
CMD ["conda", "run", "--no-capture-output", "-n", "myenv", "bash", "-c", "unzip -o model.zip && gunicorn app:app --bind 0.0.0.0:8080"]
