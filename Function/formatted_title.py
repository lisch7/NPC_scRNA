import datetime


def print_title_with_time(title):
    current_time = datetime.datetime.now()
    formatted_current_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    time_title = f"Current Time: {formatted_current_time}"
    
    total_length = 80
    equal_signs_title = "=" * ((total_length - len(title) - 2) // 2)
    title_line = f"{equal_signs_title} {title} {equal_signs_title}"

    if len(title_line) < total_length:
        title_line += "="

    equal_signs_time = "=" * ((total_length - len(time_title) - 2) // 2)
    time_line = f"{equal_signs_time} {time_title} {equal_signs_time}"

    if len(time_line) < total_length:
        time_line += "="

    formatted_title = f"\n{title_line}\n{time_line}"
    
    print(formatted_title)
    
def print_title_without_time(title):

    total_length = 80
    equal_signs = "=" * ((total_length - len(title) - 2) // 2)
    title_line = f"{equal_signs} {title} {equal_signs}"
    
    if len(title_line) < total_length:
        title_line += "="

    formatted_title = f"\n{title_line}\n"
    
    print(formatted_title)