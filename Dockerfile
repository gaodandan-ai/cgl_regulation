# Use a lightweight official Python image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /app

# Copy dependency definition
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application files
COPY . .

# Expose port 8000 (default port used by the server)
EXPOSE 8000

# Set environment variables:
# - PORT: to allow online hosting platforms (Render, Heroku, etc.) to set their custom port
# - HEADLESS: set to true to prevent run_server.py from trying to open a browser in a headless cloud server
ENV PORT=8000
ENV HEADLESS=true

# Command to run the application
CMD ["python", "run_server.py"]
