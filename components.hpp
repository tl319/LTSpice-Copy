struct node
{

};

class component
{
    protected:
        node a;
        node b;
};

class resistor: public component
{
    private:
        double resistance;
};

class capacitor: public component
{
    private:
        double capacitance;
};

class inductor: public component
{
    private:
        double inductance;
};

class voltage_source: public component
{
    private:
        double voltage;
};

class current_source: public component
{
    private:
        double current;
};

class AC_voltage_source: public component
{
    private:
        double voltage;
        double amplitude;
        double frequency;
        double phase;
};