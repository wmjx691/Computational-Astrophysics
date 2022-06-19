from PIL import Image     # install Pillow first

class Character:
    def __init__(self, name, img):
        """
        Inputs: name [string]  the name of the character
                img  [string]  the image file name of the character
        """
        self.name = name
        self.img  = img
        return

    def get_name(self):
        """
        return the name of the character
        """
        return self.name

    def setup_image(self, size):
        """
        setup the image file
        """
        image = Image.open(self.img)
        image = image.transpose(Image.ROTATE_180)
        image = image.rotate(180)
        image = image.resize((size,size))   
        self.image = image
        return

class AngryBird(Character):
    pass

class HappyPig(Character):
    pass

if __name__=='__main__':

    bird = AngryBird(name="bird",img='./images/bird.png')

    name = bird.get_name()
    print(name)

