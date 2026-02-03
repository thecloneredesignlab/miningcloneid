from PIL import Image

def combine_images_horizontally(img1_path, img2_path, output_path):
    img1 = Image.open(img1_path)
    img2 = Image.open(img2_path)

    # Create new image
    combined_width = img1.width
    combined_height = img1.height + img2.height
    combined_img = Image.new("RGB", (combined_width, combined_height))

    # Paste images
    combined_img.paste(img1, (0, 0))
    combined_img.paste(img2, (0, img1.height))

    combined_img.save(output_path)

# Example usage
combine_images_horizontally("myplot.png", "myplot1.png", "combined.png")