class Progress:

    def __init__(self, progress_bar, message):
        self.progress_bar = progress_bar
        self.message = message
        self.step = 0

    def ping(self):
        self.step += 1
        self.progress_bar.set(self.step, message=self.message)
