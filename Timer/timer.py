import tkinter as tk
import time
import csv
from datetime import datetime, timedelta

class WeeklyTimerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Weekly Work Timer")

        self.is_running = False
        self.start_time = None
        self.elapsed_time = 0

        # Create UI elements
        self.time_label = tk.Label(root, text="00:00:00", font=("Helvetica", 48), padx=20, pady=20)
        self.time_label.pack()

        button_frame = tk.Frame(root)
        button_frame.pack(pady=10)

        self.start_button = tk.Button(button_frame, text="Start", command=self.start_timer, width=12)
        self.start_button.grid(row=0, column=0, padx=5)

        self.pause_button = tk.Button(button_frame, text="Pause", command=self.pause_timer, width=12)
        self.pause_button.grid(row=0, column=1, padx=5)

        self.reset_button = tk.Button(button_frame, text="Reset", command=self.reset_timer, width=12)
        self.reset_button.grid(row=0, column=2, padx=5)

        self.save_button = tk.Button(button_frame, text="Save Record", command=self.save_record, width=12)
        self.save_button.grid(row=1, column=0, padx=5, pady=5)

        self.show_weekly_button = tk.Button(button_frame, text="Show Weekly Stats", command=self.show_weekly_stats, width=12)
        self.show_weekly_button.grid(row=1, column=1, padx=5, pady=5)

        self.show_total_time_button = tk.Button(button_frame, text="Show Total Time", command=self.show_total_time, width=12)
        self.show_total_time_button.grid(row=1, column=2, padx=5, pady=5)

        self.quit_button = tk.Button(button_frame, text="Quit", command=root.quit, width=12)
        self.quit_button.grid(row=1, column=3, padx=5, pady=5)

        self.update_timer()

    def update_timer(self):
        if self.is_running:
            now = time.time()
            self.elapsed_time = now - self.start_time
            self.display_time(self.elapsed_time)
        self.root.after(1000, self.update_timer)

    def display_time(self, elapsed_time):
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = int(elapsed_time % 60)
        time_string = f"{hours:02}:{minutes:02}:{seconds:02}"
        self.time_label.config(text=time_string)

    def start_timer(self):
        if not self.is_running:
            self.start_time = time.time() - self.elapsed_time
            self.is_running = True

    def pause_timer(self):
        if self.is_running:
            self.is_running = False

    def reset_timer(self):
        self.is_running = False
        self.elapsed_time = 0
        self.display_time(self.elapsed_time)

    def save_record(self):
        if not self.start_time:
            return
        
        end_time = time.time()
        total_time = end_time - self.start_time if self.is_running else self.elapsed_time
        date_str = datetime.now().strftime('%Y-%m-%d')
        total_time_str = f"{int(total_time // 3600):02}:{int((total_time % 3600) // 60):02}:{int(total_time % 60):02}"

        record = [date_str, total_time_str]
        self.write_to_csv(record)

    def write_to_csv(self, record):
        file_exists = False
        try:
            with open('work_records.csv', 'r') as file:
                file_exists = True
        except FileNotFoundError:
            pass
        
        with open('work_records.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            if not file_exists:
                writer.writerow(["Date", "Work Duration"])
            writer.writerow(record)

    def show_weekly_stats(self):
        week_stats = self.calculate_weekly_stats()
        stats_message = "Weekly Work Time Statistics:\n"
        for date, duration in week_stats.items():
            stats_message += f"{date}: {duration}\n"
        
        self.show_message("Weekly Statistics", stats_message)

    def calculate_weekly_stats(self):
        week_start = datetime.now() - timedelta(days=datetime.now().weekday())
        week_start_str = week_start.strftime('%Y-%m-%d')
        
        week_stats = {}
        try:
            with open('work_records.csv', 'r') as file:
                reader = csv.reader(file)
                next(reader)  # Skip header
                for row in reader:
                    date, duration = row
                    if date >= week_start_str:
                        if date not in week_stats:
                            week_stats[date] = "00:00:00"
                        week_stats[date] = self.add_times(week_stats[date], duration)
        except FileNotFoundError:
            pass

        return week_stats

    def add_times(self, time1, time2):
        h1, m1, s1 = map(int, time1.split(':'))
        h2, m2, s2 = map(int, time2.split(':'))
        
        total_seconds = (h1 * 3600 + m1 * 60 + s1) + (h2 * 3600 + m2 * 60 + s2)
        hours = total_seconds // 3600
        minutes = (total_seconds % 3600) // 60
        seconds = total_seconds % 60
        
        return f"{hours:02}:{minutes:02}:{seconds:02}"

    def show_total_time(self):
        total_time = self.calculate_total_time()
        self.show_message("Total Time This Week", f"Total Work Time This Week: {total_time}")

    def calculate_total_time(self):
        week_start = datetime.now() - timedelta(days=datetime.now().weekday())
        week_start_str = week_start.strftime('%Y-%m-%d')
        
        total_time = "00:00:00"
        try:
            with open('work_records.csv', 'r') as file:
                reader = csv.reader(file)
                next(reader)  # Skip header
                for row in reader:
                    date, duration = row
                    if date >= week_start_str:
                        total_time = self.add_times(total_time, duration)
        except FileNotFoundError:
            pass

        return total_time

    def show_message(self, title, message):
        message_box = tk.Toplevel(self.root)
        message_box.title(title)
        message_box.geometry("300x200")
        message_label = tk.Label(message_box, text=message, padx=10, pady=10)
        message_label.pack()
        close_button = tk.Button(message_box, text="Close", command=message_box.destroy)
        close_button.pack(pady=5)

# Create the main Tkinter window
root = tk.Tk()
app = WeeklyTimerApp(root)
root.mainloop()
