IMAGE_NAME="coloc-flow"
GITHUB_USERNAME="haglunda"  # Lowercase username for ghcr.io
IMAGE_TAG="latest"
GHCR_IMAGE="ghcr.io/${GITHUB_USERNAME}/${IMAGE_NAME}:${IMAGE_TAG}"

echo "Building Docker image: ${IMAGE_NAME}:${IMAGE_TAG}"
docker build -t ${IMAGE_NAME}:${IMAGE_TAG} .

echo "Tagging image for GitHub Container Registry: ${GHCR_IMAGE}"
docker tag ${IMAGE_NAME}:${IMAGE_TAG} ${GHCR_IMAGE}

echo "Please login to GitHub Container Registry"
echo "You can create a personal access token at https://github.com/settings/tokens"
echo "Token needs 'write:packages' permission"
echo "Running: docker login ghcr.io -u ${GITHUB_USERNAME}"
docker login ghcr.io -u ${GITHUB_USERNAME}

echo "Pushing image to GitHub Container Registry: ${GHCR_IMAGE}"
docker push ${GHCR_IMAGE}

echo "Image successfully pushed to GitHub Container Registry"
echo "You can pull it using: docker pull ${GHCR_IMAGE}"